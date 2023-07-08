function movavgforWG
%% set par
setup_parameters;
gsta.st=cell2mat({par.gsta(:).st}');
gsta.stn=cell2mat({par.gsta(:).stn}'); %#ok<*STRNU>
indir = './WG/';
outdir = './WGsmvrange/'; %#ok<*NASGU> 
evid = dir([indir,'20*']);
wgfileType = 'wga.R.*.nz00pct.txt';
load('./VphaseMod.mat');
parmat='smoothpar.mat';
fs = '/';
is_parfor=1;
n_parpool=8;
is_overwrite=0;
is_figure=0;
vlim = 1;% km/s, upper and lower limits of velocity in each periods
pdb=[18 50];
smrv=[0.60 0.3 0.60]; %%%smooth range for velocity
smra=[0.8  0.4 0.8];%%%smooth range for anositropy
periods=par.periods;
smrvs(periods<pdb(1))=smrv(1);
smrvs(periods>=pdb(1)&periods<=pdb(2))=smrv(2);
smrvs(periods>pdb(2))=smrv(3);
smras(periods<pdb(1))=smra(1);
smras(periods>=pdb(1)&periods<=pdb(2))=smra(2);
smras(periods>=pdb(2))=smra(3);
knst=5; % minimum number of support stations
grdz=par.gridsize;
save(parmat)

%% loop for events
if is_parfor==1
    %     parpool('local',n_parpool);
    parfor kevi = 1:length(evid) %loops for events
        smoothWG(kevi,parmat)
    end
    %     delete(gcp)
else
    for kevi = 165:length(evid) %loops for events
        smoothWG(kevi,parmat)
    end
end
delete(parmat)



%%
function smoothWG(kevi,parmat)
load(parmat)
evi = evid(kevi).name;
par.evi=evi;
wgid = dir([indir,evi,fs,wgfileType]);
for wi=1:length(wgid)    %loops for files
    disp([['kevi: ' num2str(kevi)] '  wi:' num2str(wi)])
    cntrprd = strsplit(wgid(wi).name,'.');
    cntrprd = cell2mat(cntrprd(3));
    cntrprd = str2double(cntrprd(2:4)); 
    if isempty(find(periods==cntrprd, 1))
        continue
    end
    smradv=deg2km(smrvs(periods==cntrprd));
    smradaniso=deg2km(smras(periods==cntrprd));    
    cntrV=VphaseMod(VphaseMod(:,2)==cntrprd,1); %#ok<*NODEF>
    wginpth = [indir,evi,fs,wgid(wi).name];
    outdir1 = [outdir fs evi fs];
    flo=[outdir1  wgid(wi).name num2str2(smradv,4,1) 'km.mat'];
    flot=[outdir1  wgid(wi).name num2str2(smradv,4,1) 'km.txt'];
    if is_overwrite==0&&exist(flo,'file')
        continue
    end
    %%
    [wga,weight]=wgaread(wginpth,cntrV,vlim);%
    if length(wga.st)<10||weight<0.5
        continue
    end
    wgai=VcrAmp(wga,grdz,cntrprd);
    prd.cntrprd = cntrprd;
    prd.vavg=median(wgai.v);
    kv=prd.vavg+[-vlim vlim];
    %%
    wgai=datselect(wgai,knst,kv);

    if length(wgai.st)<10
        continue
    end
    %%
    wga1=smoothavg(wgai,smradv,smradaniso);
    if isempty(wga1)
        continue
    end
    if is_figure==1&&is_parfor==0
        figure(1)
        plot_map2(wga1)
        hold off    
        figure(2)
        plot_map(wga1)
        hold off
    end

    if ~exist(outdir1,'dir')
        mkdir(outdir1)
    end
%     OutPutVelAzm(flot,wga1);
    wga=wga1;
    save(flo,'wga')
end


%
%---Read the WGA files---
%
function [wga,weight]=wgaread(fl,cntrV,vlim)

[stn,stla,stlo,stel,v,dv,az,da,rd,dr,gs,dg,vg,ao,ns,envu0,amp,Time,evla,evlo,evdpth,Ax,Ay,Bx,By,SNR]= ...
    textread(fl,'%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f'); %#ok<DTXTRD>
SNRwi=log10(SNR);
v1=max([2,cntrV-vlim]);
v2=min([5,cntrV+vlim]);
bednb=find(v>v2 | v<v1);
bednb2=find(v>5 | v<2);
weight=(1-((length(bednb)/length(v))+2*(length(bednb2)/length(v)))/3);
dltIDX=bednb;
stn(dltIDX) = [];
stla(dltIDX) = [];
stlo(dltIDX) = [];
stel(dltIDX) = [];
v(dltIDX) = [];
dv(dltIDX) = [];
az(dltIDX) = [];
da(dltIDX) = [];
rd(dltIDX) = [];
dr(dltIDX) = [];
gs(dltIDX) = [];
dg(dltIDX) = [];
vg(dltIDX) = [];
ao(dltIDX) = [];
ns(dltIDX) = [];
amp(dltIDX) = [];
evla(dltIDX) = [];
evlo(dltIDX) = [];
evdpth(dltIDX) = [];
Ax(dltIDX) = [];
Ay(dltIDX) = [];
Bx(dltIDX) = [];
By(dltIDX) = [];
envu0(dltIDX) = [];
Time(dltIDX) = [];
SNRwi(dltIDX) = [];
SNR(dltIDX)=[];

wga.st=[stla stlo stel];
wga.v=v;
wga.dv=abs(dv);
wga.az=az;
wga.da=abs(da);
wga.rd=rd;
wga.dr=abs(dr);
wga.gs=gs;
wga.dg=abs(dg);
wga.vg=vg;
wga.ao=ao;
wga.nst=ns;
wga.amp=amp;
wga.stn=stn;
wga.evla=evla;
wga.evlo=evlo;
wga.evdpth=evdpth;
wga.Ax = Ax;
wga.Ay = Ay;
wga.Bx = Bx;
wga.By = By;
wga.envu0 = envu0;
wga.Time = Time;
wga.weight=weight.*SNRwi;
wga.SNR=SNR;




%%-------------------------------
%% data selection
%%-------------------------------
function wga=datselect(wga,knst,kv)
idx=getIndx(wga,knst,kv);
wga.st=wga.st(idx,1:3);
wga.v=wga.v(idx);
wga.dv=wga.dv(idx);
wga.az=wga.az(idx);
wga.da=wga.da(idx);
wga.rd=wga.rd(idx);
wga.dr=wga.dr(idx);
wga.gs=wga.gs(idx);
wga.dg=wga.dg(idx);
wga.vg=wga.vg(idx);
wga.ao=wga.ao(idx);
wga.nst=wga.nst(idx);
wga.envu0=wga.envu0(idx);
wga.Time=wga.Time(idx);
wga.amp=wga.amp(idx);
wga.stn=wga.stn(idx);
wga.evla=wga.evla(idx);
wga.evlo=wga.evlo(idx);
wga.evdpth=wga.evdpth(idx);
wga.Ax = wga.Ax(idx);
wga.Ay = wga.Ay(idx);
wga.Bx = wga.Bx(idx);
wga.By = wga.By(idx);
wga.weight = wga.weight(idx);
wga.SNR = wga.SNR(idx);
wga.dBdx=wga.dBdx(idx);
wga.dBdy=wga.dBdy(idx);
wga.dAdy=wga.dAdy(idx);
wga.dAdx=wga.dAdx(idx);
wga.AB2=wga.AB2(idx);
wga.vcr=wga.vcr(idx);
wga.dvcr=wga.dvcr(idx);


function idx=getIndx(wga,knst,kv)
idx1=wga.v>kv(1) & wga.v<kv(2) & wga.nst>=knst;
% V95=wga.Ax;
% V1=median(V95)-1.999*std(V95);
% V2=median(V95)+1.999*std(V95);
% idx3=V95>=V1 & V95<=V2 & wga.nst;
% 
% V95=wga.Ay;
% V1=median(V95)-1.9999*std(V95);
% V2=median(V95)+1.9999*std(V95);
% idx4=V95>=V1 & V95<=V2 & wga.nst;

V95=wga.v;
V1=median(V95)-1.9995*std(V95);
V2=median(V95)+1.9995*std(V95);
idx5=V95>=V1 & V95<=V2 & wga.nst;
% 
% V95=wga.Bx;
% V1=median(V95)-1.9995*std(V95);
% V2=median(V95)+1.9995*std(V95);
% idx6=V95>=V1 & V95<=V2 & wga.nst;
% 
% V95=wga.By;
% V1=median(V95)-1.9995*std(V95);
% V2=median(V95)+1.9995*std(V95);
% idx7=V95>=V1 & V95<=V2 & wga.nst;

% V95=wga.envu0;
% V1=median(V95)-1.9999*std(V95);
% V2=median(V95)+1.9999*std(V95);
% idx8=V95>=V1 & V95<=V2 & wga.nst;
% 
% V95=wga.vcr;
% V1=median(V95)-1.9995*std(V95);
% V2=median(V95)+1.9995*std(V95);
% idx9=V95>=V1 & V95<=V2 & wga.nst;
% 
% V95=wga.dBdx;
% V1=median(V95)-1.9995*std(V95);
% V2=median(V95)+1.9995*std(V95);
% idx10=V95>=V1 & V95<=V2 & wga.nst;
% 
% V95=wga.dBdy;
% V1=median(V95)-1.9995*std(V95);
% V2=median(V95)+1.9995*std(V95);
% idx11=V95>=V1 & V95<=V2 & wga.nst;


idx=idx1 & idx5;

%%
%%average around a station
%%
function bk=smoothavg(wga,dsmax,dsmax2)
if isempty(dsmax)
    bk=wga;
    return
end
st=wga.st;
[nst,~]=size(st);
bk=[];
xo=wga.st(:,2);
ya=wga.st(:,1);
for is=1:nst
    stn=cell2mat(wga.stn(is,:));
    yi=st(is,1);
    xi=st(is,2);
    dst=distanceX(yi,xi,st(:,1),st(:,2))*111.1949;    
    idx=dst<=dsmax;
    idx2=dst<=dsmax2;
    bk.st(is,1:3)=st(is,1:3);
    bk.ao(is,1)=wga.ao(is,1);
    bk.nst(is,1)=wga.nst(is,1);
    bk.envu0(is,1)=wga.envu0(is,1);
    bk.Time(is,1)=wga.Time(is,1);
    bk.amp(is,1)=wga.amp(is,1);
    bk.stn(is,1)={stn};
    bk.evla(is,1)=wga.evla(is,1);
    bk.evdpth(is,1)=wga.evdpth(is,1);
    bk.evlo(is,1)=wga.evlo(is,1);
    bk.weight(is,1)=wga.weight(is,1);
    bk.SNR(is,1)=wga.SNR(is,1);
    bk.v(is,1)=median(wga.v(idx));
    bk.va(is,1)=median(wga.v(idx2));
    bk.dv(is,1)=median(wga.dv(idx));
    bk.az(is,1)=median(wga.az(idx));
    bk.da(is,1)=median(wga.da(idx));
    bk.rd(is,1)=median(wga.rd(idx));
    bk.dr(is,1)=median(wga.dr(idx));
    bk.gs(is,1)=median(wga.gs(idx));
    bk.dg(is,1)=median(wga.dg(idx));
    bk.Ax(is,1) = median(wga.Ax(idx));
    bk.Ay(is,1) = median(wga.Ay(idx));
    bk.Bx(is,1) = median(wga.Bx(idx));
    bk.By(is,1) = median(wga.By(idx));
    bk.dBdx(is,1)=median(wga.dBdx(idx));
    bk.dBdy(is,1)=median(wga.dBdy(idx));
    bk.dAdy(is,1)=median(wga.dAdy(idx));
    bk.dAdx(is,1)=median(wga.dAdx(idx));
    bk.AB2(is,1)=median(wga.AB2(idx));
    bk.vcr(is,1)=median(wga.vcr(idx));
    bk.vcra(is,1)=median(wga.vcr(idx2));    
    bk.dvcr(is,1)=median(wga.dvcr(idx));
end


%%
%%average around a station
%%

%-----------------------------------
%------Output results-------------
%-----------------------------------
function OutPutVelAzm(fout,wga)
nsts=length(wga.v);
if isempty(nsts)
    return
end
fido=fopen(fout,'w');
is=1:nsts;
st=wga.st(is,1:3);
st=st';
v0=wga.v(is);
va=wga.va(is);
dv=wga.dv(is);
a0=wga.az(is);
da=wga.da(is);
r0=wga.rd(is);
dr=wga.dr(is);
g0=wga.gs(is);
dg=wga.dg(is);
ao=wga.ao(is);
nst=round(wga.nst(is));
amp=wga.amp(is);
evla=wga.evla(is);
evlo=wga.evlo(is);
evdpth=wga.evdpth(is);
Ax = wga.Ax(is);
Ay = wga.Ay(is);
Bx = wga.Bx(is);
By = wga.By(is);
envu0=wga.envu0(is);
Time=wga.Time(is);
wgieht=wga.weight(is);
SNR=wga.SNR(is);
fprintf(fido,'%7.3f %8.3f %7.2f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %8.3f %2g %12.3f %12.3f %12.3f %7.3f %7.3f %7.3f %12.8f %12.8f %12.8f %12.8f %7.3f %7.3f\n', ...
    [st(1,:);st(2,:);st(3,:);v0;va;dv;a0;da;r0;dr;g0;dg;ao;nst;amp;envu0;Time;evla;evlo;evdpth;Ax;Ay;Bx;By;wgieht;SNR]);
fclose(fido);


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