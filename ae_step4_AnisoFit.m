%% Anisotropy fitting for Rayleigh wave

function ae_step4_AnisoFit
%% set up parameters
par=aa_WGM_parameters;
fs='/';

is_overwrite = 1;       % =1, overwrite results
is_parfor = 1;          % =1,multicore computing; =0,single core computing
inroot= './aniso/station_fortest/';
outpath0 = './aniso/stationAniso/'; 
outpath1 = './aniso/';
ststr='bk*.mat';

azmbin = 20;    % azimuth intervel
MaxAzm = 360;   % maxmin azimuth
is_onephase = 0;        % =1,applying ones phase term to fit anisotry
vlim = 0.6;     % km/s,  % outlier: abs(v-median(v))>vlim
is_CIR = 0; %=1,remove the outlier,[v-viso-a*cos(2*theta)-b*sin(2*theta))>0.5*vlim]
velocityType='dym';    %fitting the anisotropy using dynamic velocity (='dym')or structure veilocity(= 'strc') 
is_bootstrap = [0 1000 0.9]; % is_bootstrap(1)=1, apply bootstrap method

prdb=[15 60];
smrv=[0.2 0.2 0.2];   %deg,smoothing radius
smra=[0.35 0.3 0.35];   %deg,smoothing radius
%% ！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！



if strcmp(velocityType,'strc')
    outmatf=[outpath1 fs 'Aniso' num2str(MaxAzm) 'azmbin' num2str(azmbin) 'onephs' num2str(is_onephase) '_cr_allst.mat'];
    outmatf2=[outpath1 fs 'Aniso' num2str(MaxAzm) 'azmbin' num2str(azmbin) 'onephs' num2str(is_onephase) '_cr_smr_allst.mat'];
elseif strcmp(velocityType,'dym')
    outmatf=[outpath1 fs 'Aniso' num2str(MaxAzm) 'azmbin' num2str(azmbin) 'onephs' num2str(is_onephase) '_allst.mat'];
    outmatf2=[outpath1 fs 'Aniso' num2str(MaxAzm) 'azmbin' num2str(azmbin) 'onephs' num2str(is_onephase) '_smr_allst.mat'];
end
steps=[1,1,1]; % flag of fitting steps 
periods=par.periods; %#ok<*NASGU> %periods
fprd=periods;  %plot
pbeg=min(periods)-1;
pend=max(periods)+1;
pdb=[pbeg prdb pend]; %period boundary for smoothing radius

smrvs(length(periods))=0;
smras(length(periods))=0;
for i=1:length(pdb)-1
    id1=periods>pdb(i)&periods<=pdb(i+1)+1;
    smrvs(id1)=smrv(i);
    smras(id1)=smra(i);
end
%---------------------------------------------
bootNum = is_bootstrap(2);%
if ~exist(outpath1,'dir') 
    mkdir(outpath1)
end
if ~exist(outpath0,'dir') 
    mkdir(outpath0) 
end
sts = dir([inroot fs ststr]);
para_aniso_parfor='para_aniso_parfor.mat';
save(para_aniso_parfor);



% sts=sts(ist);
%% -----------------step 1: fit aniso
 if is_parfor==1 && steps(1)==1
    parfor ist=1:length(sts)
            parforFitAniso(ist,para_aniso_parfor)        
    end
 elseif is_parfor==0 && steps(1)==1
    for ist=1:length(sts)
        parforFitAniso(ist,para_aniso_parfor)
    end
end


%% -------------step 2: merger results
if strcmp(velocityType,'strc')
    anisoDirs = dir([outpath0 fs '/Aniso' num2str(MaxAzm) '_cr_bk_*_*.mat']);
elseif strcmp(velocityType,'dym')
    anisoDirs = dir([outpath0 fs '/Aniso' num2str(MaxAzm) '_bk_*_*.mat']);
end
if steps(2)==1
    Aniso=mergeAniso(anisoDirs,fs,outpath1,outpath0,MaxAzm,outmatf);
end
%% ---------------step 3: smooth result
load(outmatf) %#ok<*LOAD> 
Anisoin=Aniso;
save(para_aniso_parfor,'smras','smrvs','outpath1','MaxAzm','periods','Anisoin','is_overwrite');
if is_parfor==1 && steps(3)==1
    parfor ipd=1:length(periods)        
        smoothaniso(ipd,para_aniso_parfor);
    end
elseif is_parfor==0 && steps(3)==1
    for ipd=1:length(periods)        
        smoothaniso(ipd,para_aniso_parfor);
    end
end

% -----------merger smooth result
outmats=dir([outpath1 fs 'Aniso' num2str(MaxAzm) 'p*allst_smr.mat']);
for ipd=1:length(outmats)
    Anisoi=load([outpath1 fs outmats(ipd).name]);
    Anisoi=Anisoi.Aniso;
    pdx=periods==Anisoi.period;
    Aniso.viso(:,pdx) = Anisoi.viso;
    Aniso.Fai(:,pdx) = Anisoi.Fai;
    Aniso.M(:,pdx) = Anisoi.M;
    Aniso.a(:,pdx) = Anisoi.a;
    Aniso.b(:,pdx) = Anisoi.b;
    Aniso.stdv(:,pdx) = Anisoi.stdv;
    Aniso.stdFai(:,pdx) = Anisoi.stdFai;
    Aniso.stdM(:,pdx) = Anisoi.stdM;
    Aniso.stdA(:,pdx) = Anisoi.stdA;
    Aniso.stdB(:,pdx) = Anisoi.stdB;
end

save(outmatf2,'Aniso')
%% ---------------step 4: plot fitting results
    try
        load(outmatf2)
    catch
        load(outmatf)
    end
    x=Aniso.st(:,2);
    y=Aniso.st(:,1);
    xg=unique(x);
    yg=unique(y);
    [xg,yg] = meshgrid(xg,yg);
    fi=0;
    for i=fprd
        fi=fi+1;
        periods=unique(Aniso.period);
        ids=find(periods==i);
        Visos=Aniso.viso(:,ids);
        FPD=Aniso.Fai(:,ids);
        MOA=Aniso.M(:,ids)/100;
        MOA(MOA>5)=5;
        
        v1=median(Visos(~isnan(Visos)))-4*std(Visos(~isnan(Visos)));
        v2=median(Visos(~isnan(Visos)))+4*std(Visos(~isnan(Visos)));
        deleID=Visos>v2|Visos<v1;
        
        
        Visos(deleID)=nan;
        F(deleID)=nan;
        Visos(deleID)=nan;
        
        a=MOA.*sin(FPD/180*pi);
        b=MOA.*cos(FPD/180*pi);
        [n,m]=size(xg);
        for isg=1:n
            for jsg=1:m
                idx=find(x==xg(isg,jsg) & y==yg(isg,jsg));
                if isempty(idx)
                    vGrd(isg,jsg)=nan;
                    ag(isg,jsg)=nan;
                    bg(isg,jsg)=nan;
                else
                    vGrd(isg,jsg)=Visos(idx);
                    ag(isg,jsg)=a(idx);
                    bg(isg,jsg)=b(idx);
                end
                
                
            end
        end
        
        figure      
        VinGrd = griddata(y,x,Visos,yg,xg);
        [~,c]=contourf(xg,yg,VinGrd,200);
        set(c,'edgecolor','none');
        % c.LineColor='w';
        colormap(flipud(jet))
        hold on
        % caxis([-4000 4000]);
        shading interp
        colorbar
        
        hold on
        h1=quiver(xg,yg,ag,bg,'k');
        h1.ShowArrowHead = 'off';
        h2=quiver(xg,yg,-ag,-bg,'k');
        h2.ShowArrowHead = 'off';
        axis equal
        title(['T = ' num2str(i) 's'])
        hold off
    end




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -----------------Loop for periods-----------------------------------
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function parforFitAniso(kist,para_aniso_parfor)
load(para_aniso_parfor)
temp = load([inroot fs sts(kist).name]);
fnm=cell2mat(fieldnames(temp));
sti = getfield(temp,fnm); %#ok<GFLD>
disp(['ist=' num2str(kist)])

Anisoi.st = [sti.st(1),sti.st(2)];

if strcmp(velocityType,'strc')
    V=[sti.vcr];
    Va=[sti.vcra];
    outfile = [outpath0 'Aniso' num2str(MaxAzm) '_cr_' sts(kist).name];
elseif strcmp(velocityType,'dym')
    V=[sti.v];
    Va=[sti.va];
    outfile = [outpath0 'Aniso' num2str(MaxAzm) '_' sts(kist).name];
end
dazms=[sti.dazm];
if is_overwrite==0 && exist(outfile,'file')
    return
end

AZMO=[sti.azm];
stPrd=[sti.periods];
stPrd=stPrd(~isnan(stPrd));
stPrd=unique(stPrd);
for iprd = 1:length(periods)
    id=find(stPrd==periods(iprd));
    if isempty(id)
        Anisoi.period(iprd)=periods(iprd);
        Anisoi.viso(iprd) = nan;
        Anisoi.stdv(iprd) = nan;
        Anisoi.Fai(iprd) = nan;
        Anisoi.M(iprd) = nan;
        Anisoi.a(iprd) = nan;
        Anisoi.b(iprd) = nan;
        Anisoi.fitRsrm(iprd) = nan;
        Anisoi.stdFai(iprd) = nan;
        Anisoi.stdM(iprd) = nan;
        Anisoi.stdA(iprd) = nan;
        Anisoi.stdB(iprd) = nan;
        continue
    end

    Anisoi.period(iprd)=periods(iprd);
    v = V(:,id);
    va =  Va(:,id);
    dazm=dazms(:,id);
    azmo = AZMO(:,id);
    edv=abs(v-median(v));
    delID=isnan(v) | v==0 | edv>vlim | va==0 | abs(dazm)>30;
    v(delID)=[];
    va(delID)=[];
    azmo(delID)=[];
    azmo(azmo<0)=azmo(azmo<0)+360;
    if MaxAzm==180
        azmo(azmo>180)= azmo(azmo>180)-180;
    end
    
    if length(v)<5&&~isempty(v)
        Anisoi.period(iprd)=periods(iprd);
        Anisoi.viso(iprd) = median(v);
        Anisoi.stdv(iprd) = nan;
        Anisoi.Fai(iprd) = nan;
        Anisoi.M(iprd) = nan;
        Anisoi.a(iprd) = nan;
        Anisoi.b(iprd) = nan;
        Anisoi.fitRsrm(iprd) = nan;
        Anisoi.stdFai(iprd) = nan;
        Anisoi.stdM(iprd) = nan;
        Anisoi.stdA(iprd) = nan;
        Anisoi.stdB(iprd) = nan;
        continue
    elseif isempty(v)
        Anisoi.period(iprd)=periods(iprd);
        Anisoi.viso(iprd) = nan;
        Anisoi.stdv(iprd) = nan;
        Anisoi.Fai(iprd) = nan;
        Anisoi.M(iprd) = nan;
        Anisoi.a(iprd) = nan;
        Anisoi.b(iprd) = nan;
        Anisoi.fitRsrm(iprd) = nan;
        Anisoi.stdFai(iprd) = nan;
        Anisoi.stdM(iprd) = nan;
        Anisoi.stdA(iprd) = nan;
        Anisoi.stdB(iprd) = nan;
        continue
    end
    
    % remove outlier 2
if is_CIR==1
    An1 = AnisoFit(azmbin,va,azmo,360,0);
    vfit=fun1(An1.pout,azmo/180*pi);
    dv=va-vfit;
    ind=abs(dv)>vlim*0.5;
    azmo(ind)=[];
    va(ind)=[];
    v(ind)=[];
end
    %%
    if length(v)<5&&~isempty(v)
        Anisoi.period(iprd)=periods(iprd);
        Anisoi.viso(iprd) = median(v);
        Anisoi.stdv(iprd) = nan;
        Anisoi.Fai(iprd) = nan;
        Anisoi.M(iprd) = nan;
        Anisoi.a(iprd) = nan;
        Anisoi.b(iprd) = nan;
        Anisoi.fitRsrm(iprd) = nan;
        Anisoi.stdFai(iprd) = nan;
        Anisoi.stdM(iprd) = nan;
        Anisoi.stdA(iprd) = nan;
        Anisoi.stdB(iprd) = nan;
        continue
    elseif isempty(v)
        Anisoi.period(iprd)=periods(iprd);
        Anisoi.viso(iprd) = nan;
        Anisoi.stdv(iprd) = nan;
        Anisoi.Fai(iprd) = nan;
        Anisoi.M(iprd) = nan;
        Anisoi.a(iprd) = nan;
        Anisoi.b(iprd) = nan;
        Anisoi.fitRsrm(iprd) = nan;
        Anisoi.stdFai(iprd) = nan;
        Anisoi.stdM(iprd) = nan;
        Anisoi.stdA(iprd) = nan;
        Anisoi.stdB(iprd) = nan;
        continue
    end


    bandNum = MaxAzm/azmbin;    
    ig=0;
    for ibn = 1:bandNum
        bd1 = (ibn-1)*azmbin; bd2 = ibn*azmbin;
        vID = find(azmo>=bd1 & azmo<=bd2, 1);
        if isempty(vID)
            continue
        else
            ig=ig+1;
        end
    end
    %%
    if  ig<5
        Anisoi.viso(iprd) = median(v);
        Anisoi.stdv(iprd) = nan;
        Anisoi.Fai(iprd) = nan;
        Anisoi.M(iprd) = nan;
        Anisoi.a(iprd) = nan;
        Anisoi.b(iprd) = nan;
        Anisoi.fitRsrm(iprd) = nan;
        Anisoi.stdFai(iprd) = nan;
        Anisoi.stdM(iprd) = nan;
        Anisoi.stdA(iprd) = nan;
        Anisoi.stdB(iprd) = nan;
        continue
    end
    %%
    An = AnisoFit(azmbin,va,azmo,MaxAzm,is_onephase);
    Anisoi.viso(iprd) = median(v);
    Anisoi.Fai(iprd) = An.Fai;
    Anisoi.M(iprd) = An.M;
    Anisoi.fitRsrm(iprd) = An.fitRsrm;
    Anisoi.a(iprd) = An.a;
    Anisoi.b(iprd) = An.b;
    if is_bootstrap(1)==1
        fitNum = round(length(v)*is_bootstrap(3));
        bootAniso = bootaniso(fitNum,v,va,azmo,azmbin,bootNum,MaxAzm,is_onephase);
        bootFai=[bootAniso.Fai];        
        bootM=[bootAniso(:).Fai];
        boota=[bootAniso(:).a];
        bootb=[bootAniso(:).b];        
        nanID=find(isnan(bootFai));
        if length(nanID)/length(bootFai)>0.5
            Anisoi.stdFai(iprd) = nan;
            Anisoi.stdM(iprd) = nan;
            Anisoi.stdA(iprd) = nan;
            Anisoi.stdB(iprd) = nan;
            Anisoi.stdv(iprd) = nan;
            continue
        end
        bootFai(nanID)=[];
        bootM(nanID)=[];
        boota(nanID)=[];
        bootb(nanID)=[];
        dFai=bootFai-An.Fai;
        while ~isempty(find(abs(dFai)>90,1))
            bootFai(abs(dFai)>90)=bootFai(abs(dFai)>90)-180.*dFai(abs(dFai)>90)./abs(dFai(abs(dFai)>90));
            dFai=bootFai-median(bootFai);
        end
        Anisoi.stdv(iprd)=std([bootAniso.viso]);
        Anisoi.stdFai(iprd) = std(bootFai);
        Anisoi.stdM(iprd) = std([bootAniso.M]);
        Anisoi.stdA(iprd) = std([bootAniso.a]);
        Anisoi.stdB(iprd) = std([bootAniso.b]);
    else
        
        Anisoi.stdv(iprd)=nan;
        Anisoi.stdFai(iprd) = nan;
        Anisoi.stdM(iprd) = nan;
        Anisoi.stdA(iprd) = nan;
        Anisoi.stdB(iprd) = nan;
    end
end
save(outfile,'Anisoi')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------- bootstrap estimate ---------------------
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Aniso = bootaniso(fitNum,v,va,azmo,azmbin,bootNum,MaxAzm,is_onephase)
%bootstrap method
bootID = 1:fitNum;
bootID = bootstrp(bootNum, @bootline,bootID);
% [Aniso(1:bootNum).flag]=deal(1);
for iboot = 1:bootNum
    bootVa=[];
    bootV=[];
    bootazm=[];
    ID=unique(bootID(iboot,:));
    bootVa(:,1) = va(ID);
    bootV(:,1) = v(ID);
    bootazm(:,1) = azmo(ID);
    bootdat = [bootVa,bootazm];
    bootdat=sortrows(bootdat,2);
    bootVa = bootdat(:,1);
    bootazm = bootdat(:,2);
    bandNum = MaxAzm/azmbin;
    vt2(1:bandNum) = nan;
    for ibn = 1:bandNum
        bd1 = (ibn-1)*azmbin; bd2 = ibn*azmbin;
        vID = find(bootazm>=bd1 & bootazm<=bd2);
        if isempty(vID)
            continue
        else
            vt2(ibn) = median(bootVa(vID));
        end
    end
    vt=vt2(~isnan(vt2));
    if length(vt)<5
         Aniso(iboot).viso=median(bootV);
         Aniso(iboot).Fai=nan;
         Aniso(iboot).M=nan;
         Aniso(iboot).a=nan;
         Aniso(iboot).b=nan;
         Aniso(iboot).Fai1=nan;
         Aniso(iboot).M1=nan;
         Aniso(iboot).a1=nan;
         Aniso(iboot).b1=nan;
         Aniso(iboot).fitRsrm=nan;
         Aniso(iboot).pout=nan;
         continue
    end
    Aniso(iboot) = AnisoFit(azmbin,bootVa,bootazm,MaxAzm,is_onephase);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------- Anisotropy fitting ---------------------
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Aniso = AnisoFit(azmbin,v,azm,MaxAzm,is_onephase)
%          band averege
bandNum = MaxAzm/azmbin;
k=0;
for ibn = 1:bandNum
    bd1 = (ibn-1)*azmbin;
    bd2 = ibn*azmbin;
    vID = find(azm>=bd1 & azm<=bd2);
    if isempty(vID)
        continue
    else
        k = k+1;
        v2(k) = median(v(vID)); %#ok<*AGROW>
        azmb = azm(vID);
        ddv = abs(v(vID)-v2(k));
        azmID = ddv==min(ddv);
        if length(azmID)>1
            azm2 (k) = mean(azmb(azmID));
        elseif length(azmID)==1
            azm2 (k) = azmb(azmID);
        end
    end
end
azm2 = azm2/180*pi;
viso = median(v);
if k<5
    Aniso.viso = viso;
    Aniso.Fai = 0;
    Aniso.M = 0;
    Aniso.a = 0;
    Aniso.b = 0;
    Aniso.fitRsrm = nan;
    Aniso.Fai1 = 0;
    Aniso.M1 = 0;
    Aniso.a1 = 0;
    Aniso.b1 = 0; 
    Aniso.pout=[viso,0,0,0,0];
    disp('no enouth data, skip');
    return
end

if is_onephase==0 || MaxAzm==180
    lb=[viso*0.8,-0.15,-0.15];
    ub=[viso*1.2,0.15,0.15];
    co = [viso,0,0];
    options = optimoptions('lsqcurvefit','Display','off');
    [pout,anisoRsrm] = lsqcurvefit(@fun1,co,azm2,v2,lb,ub,options);
    x=0:MaxAzm;
    y=fun1(pout,x/180*pi);
    y1=y*0+viso;
    FPD=x(y==max(y));
    FPD=FPD(1);
    MOA=(max(y)-min(y))/pout(1)*100;
    FPD(FPD>180)=FPD(FPD>180)-180;
    Aniso.viso = pout(1);
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
    lb=[viso*0.8,-0.15,-0.15,-0.15,-0.15];
    ub=[viso*1.2,0.15,0.15,0.15,0.15];
    co = [viso,0,0,0,0];
    options = optimoptions('lsqcurvefit','Display','off');
    [pout,anisoRsrm] = lsqcurvefit(@fun2,co,azm2,v2,lb,ub,options);
    x=0:360;
    y=fun1(pout(1:3),x/180*pi);
    y1=fun2(pout,x/180*pi)-y+viso;
    FPD=x(y==max(y));
    FPD=FPD(1);
    FPD(FPD>180)=FPD(FPD>180)-180;
    MOA=(max(y)-min(y))/viso*100;
    FPD(FPD>180)=FPD(FPD>180)-180;
    MOA1=(max(y1)-min(y1))/viso*100;
    FPD1=x(y1==max(y1));
    Aniso.viso = pout(1);
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
Aniso.pout=pout;
return
yiso=y*0+viso; %#ok<UNRCH>
plotAnisofit(x,y,y1,yiso,azm,v,azm2,v2,FPD,MOA)
drawnow

function y1=bootline(y)
y1=y;

function fy1=fun1(abc,fx)
fy1=abc(1)+abc(2)*cos(2*fx)+abc(3)*sin(2*fx);


function fy1=fun2(abc,fx)
fy1=abc(1)+abc(2)*cos(2*fx)+abc(3)*sin(2*fx)+abc(4)*cos(fx)+abc(5)*sin(fx);



function plotAnisofit(x,y,y1,yiso,azm,v,azm2,v2,FPD,MOA) %#ok<DEFNU>
viso=mean(y);
scatter(azm,v,'ko')
hold on
scatter(azm2*180/pi,v2,40,'ro','fill')
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


function smoothaniso(ipd,parmat)
load(parmat)
pdi=periods(ipd);
smrv=smrvs(periods==pdi);
smra=smras(periods==pdi);
if isempty(smra)||smra==0
    return
end
st=Anisoin.st;
[nst,~]=size(st);
xo=st(:,2);
ya=st(:,1);
outmat=[outpath1 'Aniso' num2str(MaxAzm) 'p' num2str(pdi) 'allst_smr.mat'];
outtxf=[outpath1 'Aniso' num2str(MaxAzm) 'p' num2str(pdi) 'allst_smr.txt'];

if is_overwrite==0 && exist(outmat,'file')
    return
end

Vip=Anisoin.viso(:,ipd);
Aip=Anisoin.a(:,ipd);
Bip=Anisoin.b(:,ipd);
MOAip=Anisoin.M(:,ipd);
dvip=Anisoin.stdv(:,ipd);
dFPDip=Anisoin.stdFai(:,ipd);
dMOAip=Anisoin.stdM(:,ipd);
dAip=Anisoin.stdA(:,ipd);
dBip=Anisoin.stdB(:,ipd);

Aniso.st=st;
Aniso.period=pdi;
% ouf=fopen(outtxf,'w');
disp(ipd)
for is=1:nst
    idxv=nrGrd(st,is,smrv,xo,ya);
    idxa=nrGrd(st,is,smra,xo,ya);
    idxa(isnan(Aip(idxa))|dFPDip(idxa)>45|dMOAip(idxa)>median(MOAip(~isnan(MOAip))))=[];
    if length(idxa)<1
        Aniso.viso(is,1)=nan;
        Aniso.Fai(is,1) = 0;
        Aniso.M(is,1) = 0;
        Aniso.a(is,1) = 0;
        Aniso.b(is,1) = 0;
        Aniso.stdv(is,1)=nan;
        Aniso.stdFai(is,1) = nan;
        Aniso.stdM(is,1) = nan;
        Aniso.stdA(is,1) = nan;
        Aniso.stdB(is,1) = nan;
        continue
    end
    MOAis=MOAip(idxa);
    difM=MOAis/mean(MOAis);
    idxa(difM>1.8|difM<0.2 | isnan(Vip(idxa)))=[];
    idxv(isnan(Vip(idxv)))=[];
    if length(idxa)<1 
        Aniso.viso(is,1)=nan;
        Aniso.Fai(is,1) = 0;
        Aniso.M(is,1) = 0;
        Aniso.a(is,1) = 0;
        Aniso.b(is,1) = 0;
        Aniso.stdv(is,1)=nan;
        Aniso.stdFai(is,1) = nan;
        Aniso.stdM(is,1) = nan;
        Aniso.stdA(is,1) = nan;
        Aniso.stdB(is,1) = nan;
        continue
    end
    vis0=Vip(idxv);
    vis=Vip(idxa);
    Ais=Aip(idxa);
    Bis=Bip(idxa);
    dvis=dvip(idxa);
    dAis=dAip(idxa);
    dBis=dBip(idxa);
    dFPDis=dFPDip(idxa);
    dMOAis=dMOAip(idxa);
    
    azm=(0:360);
    vazmm=fun1([mean(vis),mean(Ais),mean(Bis)],azm*pi/180);
    cofR=zeros(size(vis));
    for i=1:length(vis)
        vazm=fun1([vis(i),Ais(i),Bis(i)],azm*pi/180);
        cofRi=corrcoef(vazm,vazmm);
        cofR(i)=cofRi(1,2); 
    end
    
    vis_m=mean(vis0);
    Ais_m=mean(Ais(cofR>=0,:));
    Bis_m=mean(Bis(cofR>=0,:));
    
    dAis_m=mean(dAis(cofR>=0,:));    
    dv_m=mean(dvis(cofR>=0,:));
    dBis_m=mean(dBis(cofR>=0,:));
    dFPDis_m=mean(dFPDis(cofR>=0,:));
    dMOAis_m=mean(dMOAis(cofR>=0,:));
    
    
    if ~isnan(Ais_m)
        vazm=fun1([mean(vis(cofR>=0,:)),Ais_m,Bis_m],azm*pi/180);
        FPD=azm(vazm==max(vazm));
        FPD=FPD(1);
        FPD(FPD>180)=FPD(FPD>180)-180;
        MOA=(max(vazm)-min(vazm))/mean(vazm)*100;
        FPD(FPD>180)=FPD(FPD>180)-180;
    else
        MOA=nan;
        FPD=nan;
    end
    
    Aniso.viso(is,1) = vis_m;
    Aniso.Fai(is,1) = FPD;
    Aniso.M(is,1) = MOA;
    Aniso.a(is,1) = Ais_m;
    Aniso.b(is,1) = Bis_m;
    
    Aniso.stdv(is,1) = dv_m;
    Aniso.stdFai(is,1) = dFPDis_m;
    Aniso.stdM(is,1) = dMOAis_m;
    Aniso.stdA(is,1) = dAis_m;
    Aniso.stdB(is,1) = dBis_m;
%     stla=ya(is,1);
%     stlo=xo(is,1);
%     if ~isnan(vis_m)
%         fprintf(ouf,'%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %5.0f\n',...
%             stla,stlo,vis_m,dv_m,FPD,dFPDis_m,MOA,dMOAis_m,Ais_m,Bis_m,pdi);
%     end
    
end
% fclose(ouf);
save(outmat,'Aniso')





function idx=nrGrd(st,is,dsmax,xs,ys)
yi=st(is,1);
xi=st(is,2);
dx=abs(xs-xi);
dy=abs(ys-yi);
stid1=find(dx<4*dsmax & dy<4*dsmax);
xn=xs(stid1);
yn=ys(stid1);

dst=ones(size(xn))*nan;
for ist=1:length(xn)
    dst(ist)=distance(yi,xi,yn(ist),xn(ist));
end

idx=dst<=dsmax;
idx=stid1(idx);


function Aniso=mergeAniso(anisoDirs,fs,outpath1,outpath0,MaxAzm,outmatf) %#ok<*INUSL>
for ia =1:length(anisoDirs)
    temp = load([outpath0 fs anisoDirs(ia).name]);
    fnm=cell2mat(fieldnames(temp));
    Anisoi = getfield(temp,fnm); %#ok<GFLD>
    Aniso.st(ia,:)=Anisoi.st;
    Aniso.viso(ia,:) = Anisoi.viso;
    Aniso.stdv(ia,:) = Anisoi.stdv;
    Aniso.Fai(ia,:) =  Anisoi.Fai;
    Aniso.M(ia,:) = Anisoi.M;
    Aniso.a(ia,:) = Anisoi.a;
    Aniso.b(ia,:) = Anisoi.b;
    Aniso.fitRsrm(ia,:) = Anisoi.fitRsrm;
    Aniso.stdFai(ia,:) = Anisoi.stdFai;
    Aniso.stdM(ia,:) = Anisoi.stdM;
    Aniso.stdA(ia,:) = Anisoi.stdA;
    Aniso.stdB(ia,:) = Anisoi.stdB;

    stla=Anisoi.st(1);
    stlo=Anisoi.st(2);
    v=Anisoi.viso;
    dv=Anisoi.stdv;
    FPD=Anisoi.Fai;
    dfpd=Anisoi.stdFai;
    MOA=Anisoi.M;
    dmoa=Anisoi.stdM;
    a=Anisoi.a;
    b=Anisoi.b;
    period=Anisoi.period;
    nanID=isnan(MOA);
    
    v(nanID)=[];
    dv(nanID)=[];
    FPD(nanID)=[];
    dfpd(nanID)=[];
    MOA(nanID)=[];
    dmoa(nanID)=[];
    a(nanID)=[];
    b(nanID)=[];
    period(nanID)=[];
%     if ~isempty(period)
%         for j=1:length(period)
%             pdi=period(j);
%             outtxtf=[outpath1 'Aniso' num2str(MaxAzm) 'p' num2str(pdi) 'allst.txt'];
%             ouf=fopen(outtxtf,'a+');                
%             fprintf(ouf,'%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %5.0f\n',...
%                 stla,stlo,v(j),dv(j),FPD(j),dfpd(j),MOA(j),dmoa(j),a(j),b(j),pdi);
%             fclose(ouf);
%         end
%     end
end
Aniso.period = Anisoi.period;
save(outmatf,'Aniso');

