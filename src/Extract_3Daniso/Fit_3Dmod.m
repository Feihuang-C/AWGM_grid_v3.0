function Fit_3Dmod
% close all
% clear all
pthi0='../';  
bk=1:5000;
% bk=136;
smr=[]; % =[],no smooth
par.MaxAzm=180;
par.isparfor=1;
par.isfigure=0; 
par.isoverwrite=1; 
par.onephase=0;
par.depthrange=[0 200];
par.outyp='aniso';% =iso,=aniso
dp=[2 25;25 50;2 50;50 130;2 130;50 100;2 100;2 190;50 190;150 190;50 150];
% dp=[100 190];
flo=[pthi0 ['mod3d_1phase' num2str(par.onephase) '_maxazm' num2str(par.MaxAzm) '.mat']];  %%4D velocity file
%% read and fitting 3D mod
if  par.isoverwrite==1 || ~exist(flo,'file')
    mod3d=Fitmod(bk,pthi0,par);
    save(flo,'mod3d');
end

%% smooth mod
load(flo);
if ~isempty(smr)
    flo=[pthi0 ['mod3d_smr' num2str(smr) '_1phase' num2str(par.onephase) '_maxazm' num2str(par.MaxAzm) '.mat']];  %%4D velocity file
    if par.isoverwrite==1 || ~exist(flo,'file')
        mod3d=smoothmod(mod3d,smr);
        save(flo,'mod3d');
    end
end
load(flo);

%% depth:mod2txt
outmod(mod3d,pthi0,par)
%% depth layers: mod2txt
for idp=1:length(dp(:,1))  
    disp(idp)
    dp2=dp(idp,:);
    out_gmt_drange(mod3d,dp2,pthi0,par)
end

function outmod(vmod,pthi0,par)
ndp=length(vmod(1).dp);
lo=[vmod.xbk];
la=[vmod.ybk];
vs=cell2mat({vmod.viso}');
as=cell2mat({vmod.a}');
bs=cell2mat({vmod.b}');
FPDs=cell2mat({vmod.FPD}');
MOAs=cell2mat({vmod.MOA}');
depths=vmod(1).dp;
for idp=1:length(vs(1,:))
    dpth=depths(idp);
    vsi=vs(:,idp);
    FPD=FPDs(:,idp);
    MOA=MOAs(:,idp);
    
    fl=[pthi0 'vsdep_' num2str(dpth) '_km_' num2str(par.MaxAzm) '.' num2str(par.onephase) '.dat'];
    fid=fopen(fl,'w');
    for i=1:length(vsi)
        fprintf(fid,'%7.3f %7.3f %7.3f %7.3f %7.1f\n',la(i),lo(i),vsi(i),MOA(i),FPD(i));
    end;
    fclose(fid);
end

%% output to GMT format
function out_gmt_drange(mod3d,dp2,pthi0,par)
onephase=par.onephase;
MaxAzm=par.MaxAzm;
isfigure=par.isfigure;
nbk=length(mod3d);
ndp=length(mod3d(1).dp);
dp=mod3d(1).dp;
%output layer by layer
    di=dp2(1);
    sbk1=[];
    if(di<10)
        sbk1=['000' num2str(di)];
    elseif(di<100)
        sbk1=['00' num2str(di)];
    elseif(di<1000)
        sbk1=['0' num2str(di)];
    else
        sbk1= num2str(di);
    end        
    di=dp2(2);
    sbk2=[];
    if(di<10)
        sbk2=['000' num2str(di)];
    elseif(di<100)
        sbk2=['00' num2str(di)];
    elseif(di<1000)
        sbk2=['0' num2str(di)];
    else
        sbk2=num2str(di);
    end 
    
    fout=[pthi0 'vslay_dep' sbk1 '_' sbk2 'km_' num2str(MaxAzm) '.' num2str(onephase) '.dat'];

    di1=dp2(1);
    di2=dp2(2);
    idx=[];
    idx(1,:)=find(dp>=di1 & dp<=di2);
    for ibk=1:nbk
        xo(ibk)=mod3d(ibk).xbk;
        ya(ibk)=mod3d(ibk).ybk;
        vs=mod3d(ibk).vs;
        azm=mod3d(ibk).azm(1,:);
        if(length(idx)==1)
            vsi=vs(idx,:);
        else
            dpi=dp(idx);
            vsi0=vs(idx,:);
            ddpi=diff(dpi);
            ddpi=[ddpi;ddpi(length(dpi)-1)];
            [mc,nc]=size(vsi0);
            for ic=1:nc
                [t(:,ic),vsi(ic)]=v_eff(ddpi,vsi0(:,ic)); %%using equivelent velocity
            end
            for jc=1:mc
                ti=t(jc,:);
                dt(jc)=max(ti)-min(ti);
            end
        end
        
        [~,nc]=size(vsi);
        if nc<4
            fpd(ibk)=0;
            moa(ibk)=0;
            dltime(ibk)=0;
            viso(ibk)=mean(vsi);
        else
            Aniso=AnisoFit(vsi,azm,180,0,isfigure);
            fpd(ibk)=Aniso.Fai;
            moa(ibk)=Aniso.M;
            viso(ibk)=mean(vsi);
            dltime(ibk)=sum(dt);
        end
        kbk(ibk)=mod3d(ibk).ibk;
        %%find the parameters by fitting SINE(2theta) function
    end
    fid=fopen(fout,'w');
    fprintf(fid,'%7.3f %7.3f %6.3f %7.2f %7.1f %7.3f %d\n',[ya' xo' viso' moa' fpd' dltime' kbk']');
    fclose(fid);

    %%
%% equvilent velocity
%%
function [dt,vef]=v_eff(dp,vs)
dt=dp./vs;
vv=sum(vs.*vs.*dt)/sum(dt);
vef=sqrt(vv);


%%
%% fit to find the best parameters
%% AZD--azimuth in degree  VI--velocity
%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------- Anisotropy fitting ---------------------
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Aniso = AnisoFit(v,azm,MaxAzm,is_onephase,isfigure)

if MaxAzm==180;
azm(azm>360)=azm(azm>360)-360;
azm(azm<0)=azm(azm<0)+360;    
azm(azm>180)=azm(azm>180)-180;
azm(azm<0)=azm(azm<0)+180;
end


v2=v;
azm=azm/180*pi;
azm2=azm;
viso = median(v);
if is_onephase==0 || MaxAzm==180
    lb=[viso*0.8,-0.5,-0.5];
    ub=[viso*1.2,0.5,0.5];
    co = [viso,0,0];
    options = optimoptions('lsqcurvefit','Display','off');
    [pout,anisoRsrm] = lsqcurvefit(@fun1,co,azm2,v2,lb,ub,options);
    x=0:MaxAzm;
    y=fun1(pout,x/180*pi);
    y1=y*0+viso;
    FPD=x(y==max(y));
    FPD=FPD(1);
    MOA=(max(y)-min(y))/viso*100;
    FPD(FPD>180)=FPD(FPD>180)-180;
    Aniso.Fai = FPD;
    Aniso.M = MOA;
    Aniso.a = pout(2);
    Aniso.b = pout(3);
    Aniso.fitRsrm = anisoRsrm;
    Aniso.Fai1 = nan;
    Aniso.M1 = nan;
    Aniso.a1 = nan;
    Aniso.b1 = nan;
else
%     up2=[viso,0.1,0.1,0.005,0.005];
%     v=fun2(up2,azm2);
%     n=(rand(1,length(azm2))*2-1.0)*0.2;
%     v2=n+v;
    lb=[viso*0.8,-0.5,-0.5,-0.5,-0.5];
    ub=[viso*1.2,0.5,0.5,0.5,0.5];
    co = [viso,0,0,0,0];
    options = optimoptions('lsqcurvefit','Display','off');    
    [pout,anisoRsrm] = lsqcurvefit(@fun2,co,azm2,v2,lb,ub,options);
    x=0:MaxAzm;
    y=fun1(pout(1:3),x/180*pi);
    y1=fun2(pout,x/180*pi)-y+viso;
    FPD=x(y==max(y));
    FPD=FPD(1);
    FPD(FPD>180)=FPD(FPD>180)-180;
    MOA=(max(y)-min(y))/viso*100;
    FPD(FPD>180)=FPD(FPD>180)-180;
    MOA1=(max(y1)-min(y1))/viso*100;
    FPD1=x(y1==max(y1));
    Aniso.Fai = FPD;
    Aniso.M = MOA;
    Aniso.a = pout(2);
    Aniso.b = pout(3);
    Aniso.fitRsrm = anisoRsrm;
    Aniso.Fai1 = FPD1;
    Aniso.M1 = MOA1;
    Aniso.a1 = pout(4);
    Aniso.b1 = pout(5);    
end
Aniso.viso = pout(1);

if isfigure==1
yiso=y*0+viso; 
figure(887)
plotAnisofit(x,y,y1,yiso,azm,v,azm2,v2,FPD,MOA)
drawnow
end


function fy1=fun1(abc,fx)
fy1=abc(1)+abc(2)*cos(2*fx)+abc(3)*sin(2*fx);


function fy1=fun2(abc,fx)
fy1=abc(1)+abc(2)*cos(2*fx)+abc(3)*sin(2*fx)+abc(4)*cos(fx)+abc(5)*sin(fx);


function plotAnisofit(x,y,y1,yiso,azm,v,azm2,v2,FPD,MOA) 
viso=mean(y);
azm=azm*180/pi;
azm2=azm2*180/pi;
scatter(azm,v,'ko')
hold on
scatter(azm2,v2,40,'ro','fill')
plot(x,y1,'g','LineWidth',1.5);
plot(x,yiso,'m','LineWidth',1.5);
plot(x,y+y1-viso,'k','LineWidth',1.5);
plot(x,y,'b','LineWidth',1.5);

N1=viso-0.8;
k1=roundn(N1,-2) ;
N2=viso+0.8;
k2=roundn(N2,-2) ;
ylim([k1 k2]);
set(gca,'YTick',(k1:0.2:k2));
set(gca,'LineWidth',1.5,'FontSize',14)
xlabel('amz')
ylabel('v_p_h_s(Km/s)')
box on
title([FPD MOA])
hold off






%%
%% Read in all models
%%
function mod4ds=Fitmod(bk,pthi0,par)
depthrange=par.depthrange;
isparfor=par.isparfor;
onephase=par.onephase;
MaxAzm=par.MaxAzm;
isfigure=par.isfigure;
isverwrite=par.isoverwrite;
outyp=par.outyp;
if isparfor==1
    parfor ibk=1:length(bk)
        isfigure=0;
        fitmodi(pthi0,ibk,par)
    end
elseif isparfor==0
    for ibk=1:length(bk)
        fitmodi(pthi0,ibk,par)
    end
end
ibka=0;
for ibk=1:length(bk)
    disp(ibk)
    if(ibk<10)
        sbk=['0000' num2str(ibk)];
    elseif(ibk<100)
        sbk=['000' num2str(ibk)];
    elseif(ibk<1000)
        sbk=['00' num2str(ibk)];
    else
        sbk=['0' num2str(ibk)];
    end
    pthi=[pthi0 sbk '/'];
    if(exist([pthi 'mod4d.mat'],'file'))
        load([pthi 'mod4d.mat']);
        %%block idex
        ibka=ibka+1;          
        mod4ds(ibka).ibk=mod2d.ibk; %#ok<*AGROW>
        mod4ds(ibka).xbk=mod2d.xbk;
        mod4ds(ibka).ybk=mod2d.ybk;
        mod4ds(ibka).dp=mod2d.dp;
        mod4ds(ibka).azms=mod2d.azms;        
        mod4ds(ibka).viso=mod2d.viso; 
        mod4ds(ibka).FPD=mod2d.FPD;
        mod4ds(ibka).MOA=mod2d.MOA;
        mod4ds(ibka).a=mod2d.a;
        mod4ds(ibka).b=mod2d.b;
        for i=1:length(mod2d.a)
            mod4ds(ibka).vs(i,:)=fun1([mod2d.viso(i),mod2d.a(i),mod2d.b(i)],(0:10:170)/180*pi);
            mod4ds(ibka).azm(i,:)=0:10:170;
        end
    end
end

function  fitmodi(pthi0,ibk,par)
depthrange=par.depthrange;
isparfor=par.isparfor;
onephase=par.onephase;
MaxAzm=par.MaxAzm;
isfigure=par.isfigure;
isverwrite=par.isoverwrite;
outyp=par.outyp;
disp(num2str(ibk))
if(ibk<10)
    sbk=['0000' num2str(ibk)];
elseif(ibk<100)
    sbk=['000' num2str(ibk)];
elseif(ibk<1000)
    sbk=['00' num2str(ibk)];
else
    sbk=['0' num2str(ibk)];
end
pthi=[pthi0 sbk '/'];
outf=[pthi 'mod4d.mat'];
if(exist([pthi 'xyloc.txt'],'file'))
    floc=[pthi 'xyloc.txt'];
    nxy=load(floc);
    
    [dp,vs,azms]=readmod(pthi,sbk,depthrange);
    isovw=isverwrite==0 & exist(outf,'file');

    if(isempty(dp))||isovw
        return
    end
    
    if strcmp(outyp,'aniso')
        if length(azms)<5
            return
        end
        for dpi=1:length(dp)
            Aniso = AnisoFit(vs(dpi,:),azms,MaxAzm,onephase,isfigure);
            a(dpi)=Aniso.a;
            b(dpi)=Aniso.b;
            viso(dpi)=Aniso.viso;
            FPD(dpi)=Aniso.Fai;
            MOA(dpi)=Aniso.M;
        end
        
    elseif  strcmp(outyp,'iso')
        for dpi=1:length(dp)
            a(dpi)=0;
            b(dpi)=0;
            viso(dpi)=median(vs(dpi,:));
            FPD(dpi)=0;
            MOA(dpi)=0;
        end
    end
    %%block idex
    mod2d.ibk=nxy(1); %%block idex
    mod2d.xbk=nxy(2); %%block lon
    mod2d.ybk=nxy(3); %%block lat
    nd=length(dp);
    mod2d.dp=dp;
    mod2d.azms=azms;
    mod2d.FPD=FPD;
    mod2d.MOA=MOA;
    mod2d.viso=viso;
    mod2d.a=a;
    mod2d.b=b; %#ok<STRNU>
    save(outf,'mod2d')
end
%%
%%read in model inverted by CPS
%%
function [dp,vs,azms]=readmod(pthi,bk,depthrange)
ia=0;
 
fls=dir([pthi 'bk.' bk '.deg.*.dsp.modo3']);
if isempty(fls)
    dp=[];
    vs=[];
    azms=[];
    disp('no mod are readed')
    return
end

for i=1:length(fls)
    azmf=fls(i).name;
    azm=strsplit(azmf,'.');
    azmi=str2double(cell2mat(azm(4)));
    finv=[pthi '/' fls(i).name];  
    s=dir(finv);
    if s.bytes==0
        continue
    end
    ia=ia+1;
    azms(ia)=azmi;
    fid=fopen(finv,'r');
    str=fgets(fid);
    ii=strfind(str,'H(KM)');
    while isempty(ii)
        str=fgets(fid);
        ii=strfind(str,'H(KM)');
    end
    ih=0;
    dp0=0;
    while ~feof(fid) && dp0>=depthrange(1) && dp0<=depthrange(2)
        ih=ih+1;
        vv=fscanf(fid,'%g\n',[10]);
        vs(ih,ia)=vv(3);  %%VS
        dp(ih,1)=dp0+vv(1); %%dep
        dp0=dp(ih,1);
    end
    fclose(fid);
end
dp(end)=dp(end)+0.1;

