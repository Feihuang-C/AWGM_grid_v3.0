function ad_step3_event2sta
%% setup parameters
par=aa_WGM_parameters;
fs='/';
datarootv='./WG/';
EVS = dir([datarootv fs '2007*']); %#ok<*NASGU>
fstr='wga.smr.R*.p';
outpath='./aniso/stations_prd/';
outpath2='./aniso/station/'; 
is_parfor=1; %=1,parallel
is_overwrite=1;%=1,overwrite; =0,no overwrite

step=[1,1]; %choose processing steps, default value is [1,1]

if ~exist(outpath,'dir')
    mkdir(outpath)
end

periods = par.periods;
ev2stapar='ev2st_par0.mat';
save(ev2stapar);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 1: combine results of all earthquakes for each grid point
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(find(step==1, 1))
    if is_parfor==1
        parfor pdi=1:length(periods)
            ev2staparfor(pdi,ev2stapar)
        end
    else
        for pdi=1:length(periods)
            ev2staparfor(pdi,ev2stapar)
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step 2: marge results of all periods for each grid point
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(find(step==1, 1))
    gsta=par.gsta;
    yx=cell2mat({gsta.st}');
    xo=yx(:,2);
    ya=yx(:,1);
    inpath=outpath;
    outpath=outpath2;
    if ~exist(outpath,'dir')
        mkdir(outpath)
    end
    save(ev2stapar);
    if is_parfor==1
        parfor ist=1:length({gsta.stn}') 
            prd2staparfor(ist,ev2stapar)
        end
    else
        for ist=1:length({gsta.stn}')
            prd2staparfor(ist,ev2stapar)
        end
    end
end
delete(ev2stapar)





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------------for one earthquake-----------------------------------
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ev2staparfor(pdi,ev2stapar)
load(ev2stapar)
period=periods(pdi);
pdir=[outpath fs 'p' num2str2(period,3,0)];
pstr=num2str2(period,3,0);
logf=[outpath fs 'p' pstr '_loops.log'];
if ~exist(pdir,'dir')
    mkdir(pdir)
end
for evi=1:length(EVS)
%     if exist(logf,'file') 
%         [done_pi,done_ev]=textread(logf,'%f %f'); %#ok<*DTXTRD>
%         flgloop=find(done_pi==pdi&done_ev==evi, 1);
%         if ~isempty(flgloop)
%             continue
%         end
%     end
%     loopfid=fopen(logf,'a+');
%     fprintf(loopfid,'%g %g\n', pdi, evi);
%     fclose(loopfid);    
    
    disp(['T=' pstr '; EQ=' num2str(evi)])
    inpath1=[datarootv fs EVS(evi).name];
    infile=dir([inpath1 fs fstr pstr '*.mat']);
    if isempty(infile)
        continue
    end   
    infile1=[inpath1 fs infile(1).name];
    Temp=load(infile1);
    strcTemp=fieldnames(Temp);
    wga=getfield(Temp,cell2mat(strcTemp)); %#ok<*GFLD>
    x1=wga.st(:,2);
    y1=wga.st(:,1);
    for is=1:length(x1)
        xi=x1(is);
        yi=y1(is);
        stn=['bk_' num2str2(xi,6,2) '_' num2str2(yi,5,2)];
        st=wga.st(is,1:3);
        v0=wga.v(is);
        va=wga.va(is);
        vcr=wga.vcr(is);
        vcra=wga.vcra(is);       
        dcrv=wga.dvcr(is);        
        azmo=wga.azmo(is);
        dazm=wga.az(is);
        evla=wga.evla(is);
        evlo=wga.evlo(is);
        Bx = wga.Bx(is);
        By = wga.By(is);
        wgieht=wga.weight(is);
        outfile=[pdir fs 'T_' pstr,'_',stn '.txt' ];
        fidst=fopen(outfile,'a+');
        fprintf(fidst,'%7.3f %7.3f  %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %9.6f %9.6f %7.3f\n', ...
            st(1),st(2),st(3),period,v0,va,vcr,vcra,dcrv,azmo,dazm,evla,evlo,Bx,By,wgieht);
        fclose(fidst);
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -------------------for one period-------------------------------------
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  prd2staparfor(ist,ev2stapar)
disp(num2str2(ist,5,0))
load(ev2stapar)
xi=xo(ist);
yi=ya(ist);
stn=['bk_' num2str2(xi,6,2) '_' num2str2(yi,5,2)];
% if exist([outpath fs stn '.mat'],'file')&&is_overwrite==0
%     return
% end
kflg=0;
for pdi = 1:length(periods)
    period=periods(pdi);
    pstr=num2str2(period,3,0);
    pdir=[inpath fs 'p' pstr];
    fl=[pdir fs 'T_' pstr,'_',stn '.txt'];
    if ~exist(fl,'file')
        wgst.v(:,pdi) =nan;
        wgst.vcr(:,pdi) =nan;
        wgst.va(:,pdi) =nan;
        wgst.vcra(:,pdi) =nan;        
        wgst.dvcr(:,pdi) =nan;        
        wgst.weight(:,pdi)=nan;
        wgst.azm(:,pdi) =nan;
        wgst.dazm(:,pdi) =nan;
        wgst.rd(:,pdi) =nan;
        wgst.Bx(:,pdi) =nan;
        wgst.By(:,pdi) =nan;
        continue
    end
    [stla,stlo,stel,~,v,va,vcr,vcra,dvcr,azm,dazm,evla,evlo,Bx,By,weight]= ...
        textread(fl,' %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f'); %#ok<*DTXTRD>
    if pdi==1||kflg==0
        kflg=kflg+1;
        wgst.periods(pdi)=periods(pdi);
        wgst.st = [stla(1),stlo(1),stel(1)];
        i=1:length(v);
        wgst.evla(i,1)=evla;
        wgst.evlo(i,1)=evlo;
        wgst.v(i,pdi) =v;
        wgst.va(i,pdi) =va;
        wgst.vcr(i,pdi) =vcr;
        wgst.vcra(i,pdi) =vcra;
        wgst.dvcr(i,pdi) =dvcr;
        wgst.weight(i,pdi)=weight;
        wgst.azm(i,pdi) =azm;        
        wgst.dazm(i,pdi) =dazm;        
        wgst.Bx(i,pdi) =Bx;
        wgst.By(i,pdi) =By;
    else
        wgst.periods(pdi)=periods(pdi);
        a=complex(wgst.evla,wgst.evlo);
        b=complex(evla,evlo);

        [~, ia, ib] = intersect(a,b);
        [~,ic]=setdiff(b,a);
        [~,id]=setdiff(a,b);
        
        %% same ev 
        wgst.evla(ia)=evla(ib);
        wgst.evlo(ia)=evlo(ib);
        wgst.v(ia,pdi) =v(ib);
        wgst.va(ia,pdi) =va(ib);
        wgst.vcr(ia,pdi) =vcr(ib);
        wgst.vcra(ia,pdi) =vcra(ib);
        wgst.dvcr(ia,pdi) =dvcr(ib);
        wgst.weight(ia,pdi)=weight(ib);
        wgst.azm(ia,pdi) =azm(ib);
        wgst.dazm(ia,pdi) =dazm(ib);
        wgst.Bx(ia,pdi) =Bx(ib);
        wgst.By(ia,pdi) =By(ib);

        if ~isempty(ic)
            pdin1=pdi-1;
            iad=(length(wgst.v(:,1))+1):(length(wgst.v(:,1))+length(ic));
            wgst.v(iad,pdin1) =nan;
            wgst.va(iad,pdin1) =nan;
            wgst.vcr(iad,pdin1) =nan;
            wgst.vcra(iad,pdin1) =nan;
            wgst.dvcr(iad,pdin1) =nan;
            wgst.weight(iad,pdin1)=nan;
            wgst.azm(iad,pdin1) =nan;           
            wgst.dazm(iad,pdin1) =nan;           
            wgst.Bx(iad,pdin1) =nan;
            wgst.By(iad,pdin1) =nan; 
            
            wgst.evla(iad)=evla(ic);
            wgst.evlo(iad)=evlo(ic);
            wgst.v(iad,pdi) =v(ic);
            wgst.va(iad,pdi) =v(ic);  
            wgst.vcr(iad,pdi) =vcr(ic);
            wgst.vcra(iad,pdi) =vcra(ic);
            wgst.dvcr(iad,pdi) =dvcr(ic);
            wgst.weight(iad,pdi)=weight(ic);
            wgst.azm(iad,pdi) =azm(ic);           
            wgst.dazm(iad,pdi) =dazm(ic);           
            wgst.Bx(iad,pdi) =Bx(ic);
            wgst.By(iad,pdi) =By(ic);
        end
        if ~isempty(id)
            iad=id;            
            wgst.v(iad,pdi) =nan;
            wgst.va(iad,pdi) =nan;        
            wgst.vcr(iad,pdi) =nan;
            wgst.vcra(iad,pdi) =nan;
            wgst.dvcr(iad,pdi) =nan;
            wgst.weight(iad,pdi)=nan;
            wgst.azm(iad,pdi) =nan;           
            wgst.dazm(iad,pdi) =nan;           
            wgst.Bx(iad,pdi) =nan;
            wgst.By(iad,pdi) =nan; 
        end
    end    
end
v=wgst.v;
v=v(~isnan(v)|v>0);
if isempty(v)||kflg==0
    return
end

[zroID] = find(wgst.v==0);
wgst.v(zroID) =nan;
wgst.va(zroID) =nan;
wgst.vcr(zroID) =nan;
wgst.vcra(zroID) =nan;
wgst.dvcr(zroID) =nan;
wgst.weight(zroID)=nan;
wgst.azm(zroID) =nan;
wgst.dazm(zroID) =nan;
wgst.Bx(zroID) =nan;
wgst.By(zroID) =nan;
save([outpath fs stn '.mat'],'wgst');

